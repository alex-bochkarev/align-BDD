"""
UFL -- Uncapacitated Facility Location
(testing align-BDD machinery for specific applications)

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from graphviz import Digraph
import BDD as DD

import numpy as np

m = 4; n = 2

S = []

## define a toy problem
S.append([1]); S.append([1,2]); S.append([1,2]); S.append([2])
F = [1, 1]
λ = 0.5

## draw a problem diagram
dot = Digraph('G',comment="Uncapacitated Facility Location")
for i in range(n):
    dot.node(f"F{i+1}", shape='doublecircle', style='filled', color='lightgrey')

for j, S_j in enumerate(S):
    dot.node(f"C{j}", shape='circle')
    for i in S_j:
        dot.edge(f"F{i}", f"C{j}")

dot.view()

print([f"z{i}-{j+1}" for j in range(m) for i in S[j]])

print(len(np.sum(S)))

def create_covering_BDD(S,n,m):
    ## 1) create covering BDD
    C = DD.BDD(
        N=len(np.sum(S)),
        vars=[f"z{i}-{j+1}" for j in range(m) for i in S[j]]
    )

    root = C.addnode(None)

    # current nodes
    s_maybe = root
    s_no = None
    # TODO: assert all Sj != []

    for j in range(m):
        s_yes = None
        if j == m-1:
            S_j = S[j][:-1]
        else:
            S_j = S[j]

        for idx, i in enumerate(S_j):
            print(f"j={j},idx={idx}, maybenode={s_maybe.id}")
            new_maybe = C.addnode(s_maybe, "lo")
            new_yes = C.addnode(s_maybe, "hi")

            if not (s_yes is None):
                s_yes.link(new_yes, "hi")
                s_yes.link(new_yes, "lo")

            if not (s_no is None):
                if idx < len(S[j])-1:
                    new_no = C.addnode(s_no)
                else:
                    new_no = new_maybe
                    s_no.link(new_no, "hi")

                s_no.link(new_no, "lo")
                s_no = new_no

            s_yes = new_yes

            if idx == len(S[j])-1:
                # the last layer of the subnetwork
                # (related to customer j)
                s_no = new_maybe
                s_maybe = new_yes
                s_yes.link(new_yes, "hi")
                s_yes.link(new_yes, "lo")
                print(f"LAST F: yes:{s_yes.id}, no:{s_no.id}, maybe:{s_maybe.id}")
            else:
                s_maybe = new_maybe

    s_maybe.link(C.nodes[DD.NTRUE], "hi")
    s_maybe.link(C.nodes[DD.NFALSE], "lo")

    if not (s_yes is None) and (s_yes.id != s_maybe.id):
        s_yes.link(C.nodes[DD.NTRUE], "hi")
        s_yes.link(C.nodes[DD.NTRUE], "lo")

    if not (s_no is None):
        s_no.link(C.nodes[DD.NFALSE], "hi")
        s_no.link(C.nodes[DD.NFALSE], "lo")

    return C

create_covering_BDD(S, n, m).dump_gv().view()
## 2) create availability BDD
