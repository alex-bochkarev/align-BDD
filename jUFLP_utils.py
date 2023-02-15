"""Utility / helper functions for j-UFLP: used save / load / draw instances."""
import json
from jUFLP_cavemen import jUFLPEncoder


def save_inst(i1, i2, join_map, filename):
    """Saves the jUFLP instance to ``.json`` file.

    The file contains the following objects

    - ``inst1`` and ``inst2`` -- description of the two
      sub-instances, with identical structure containing
      the following records:

       - `S`: a list of neighborhood lists,
       - `f`: a list of lists of overlap costs,
       - `c`: a list of facility location costs.
       - `caves`: a parameter describing the points allocated to each "cave"
         (was used in legacy experiments)
    """
    with open(filename, "w") as fout:
        fout.write(json.dumps({
            'inst1':{
                'S':i1[0],
                'f':i1[1],
                'c':i1[2],
                'caves': i1[3]},
            'inst2':{
                'S':i2[0],
                'f':i2[1],
                'c':i2[2],
                'caves': i2[3]},
            'jmap':{int(j): join_map[j] for j in join_map}},
                              cls=jUFLPEncoder))


def load_inst(filename):
    """Loads a jUFLP instance from ``.json`` file.

    Returns:
      [[S,f,c,caves], [S2,f2,c2,caves], jmap]
    """
    with open(filename, "r") as fin:
        json_inst = fin.read()

    json_dct = json.loads(json_inst)
    return [[json_dct[f'inst{i}']['S'],
            json_dct[f'inst{i}']['f'],
            json_dct[f'inst{i}']['c'],
            json_dct[f'inst{i}']['caves']]
            for i in [1, 2]] + [{int(j1): json_dct['jmap'][j1]
                                 for j1 in json_dct['jmap']}]


def draw_jUFLP_inst(i1, i2, link, filename="tmp/jUFLP.dot"):
    """Saves an instance to a ``.dot`` file."""
    with open(filename, "w") as fout:
        fout.write("graph G {\n")

        for (inst, no, pref) in [(i1, 1, 'f'), (i2, 2, 's')]:
            S, f, c, caves = inst
            fout.write(f"    subgraph cluster_{no-1}" +
                       " {\n")  # encoding a sub-instance
            fout.write(f'        color=blue; label="sub-UFLP-{no}";')
            added = set([])
            for i in range(len(S)):
                for j in S[i]:
                    if ((i+1) != j) and not (((j, (i+1)) in added)
                                             or ((i+1, j) in added)):

                        fout.write(f"        {pref}{i+1}--{pref}{j};\n")
                        added.add(((i+1), j))

            fout.write("    };\n")  # end of sub-instance

        for j in link:
            fout.write(f"    f{j} -- s{link[j]}[color=red, style=dashed, penwidth=1];\n")
        fout.write("}")
