# Heap implementation (with remove)
# A. Bochkarev, 07-2019
# based on SO:
# https://stackoverflow.com/questions/13800947/deleting-from-python-heapq-in-ologn/13801331
from heapq import heappop, heappush, heapify
import copy as copy

## heap implementation
## available functions:
##   push   -- pushes an elem into the heap
##   pop    -- pops the smalles element
##   remove -- removes a specific element (lazy-deletion)
##   root   -- returns the smallest element (without popping)
class heap:
    def __init__(self, values=None):
        if values is None:
            self.heap = []
            self.N = 0
            self.tuple_type = None
        else:
            self.heap = copy.copy(values)
            heapify(self.heap) # in-place, linear time
            self.N = len(values)
            if type(values[0])==tuple:
                self.tuple_type = True
            else:
                self.tuple_type = False

        self.tombstones = {}

    # we have to be careful since the first element can be deleted
    def root(self):
        if len(self.heap)==0:
            return None

        while self.tombstones.get(self.heap[0],0) > 0:
            self.tombstones[self.heap[0]] -= 1
            heappop(self.heap)
            if len(self.heap)==0:
                return None

        return self.heap[0]

    def __str__(self):
        ts = copy.copy(self.tombstones)
        els = []
        if self.tuple_type:
            for i in range(len(self.heap)):
                if ts.get(self.heap[i][0],0) == 0:
                    els.append(self.heap[i])
                else:
                    ts[self.heap[i][0]] -= 1
        else:
            for i in range(len(self.heap)):
                if ts.get(self.heap[i],0) == 0:
                    els.append(self.heap[i])
                else:
                    ts[self.heap[i]] -= 1

        return "{}".format(els)

    # returns the no. of ("alive") elements
    def __len__(self):
        return self.N

    def push(self, what):
        if type(what)==tuple:
            self.tuple_type = True
        else:
            self.tuple_type = False

        heappush(self.heap, what)
        self.N += 1

    def remove(self, what):
        if self.tuple_type:
            self.tombstones[what] = self.tombstones.get(what,0)+1 # mark the key for deletion
        else:
            self.tombstones[what] = self.tombstones.get(what,0)+1 # mark for deletion

        self.N -= 1
        if self.N<0:
            self.N = 0
            print("Error: the heap is empty (attempted to delete {})".format(what))

    # pops the smallest element
    # where all the magic happens
    def pop(self):
        # first, clean up all the dead bodies
        if self.tuple_type:
            while len(self.heap) > 0 and self.heap[0][0] in self.tombstones and self.tombstones[self.heap[0][0]]:
                self.tombstones[self.heap[0][0]] -= 1
                heappop(self.heap)
        else:
            while len(self.heap) > 0 and self.heap[0] in self.tombstones and self.tombstones[self.heap[0]]:
                self.tombstones[self.heap[0]] -= 1
                heappop(self.heap)

        # now, pop for real
        if len(self.heap)>0:
            self.N -= 1
            return heappop(self.heap)
        else:
            print("Error: the heap is empty")
            return None
