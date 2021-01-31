"""
UFL -- Uncapacitated Facility Location
(testing align-BDD machinery for specific applications)

---
(c) A. Bochkarev, Clemson University, 2021
abochka@clemson.edu
"""
from graphviz import Digraph

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
