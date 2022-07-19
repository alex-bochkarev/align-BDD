"""A simple illustration for the special instance type for joint UFLP."""
from UFLP_2_cav import gen_special_jUFLP
from jUFLP_utils import draw_jUFLP_inst
from UFLPOrder import UFLP_greedy_order, N2RList
from UFLP_fullDD import create_cover_DD
from jUFLP_utils import load_inst, save_inst


# if __name__ == '__main__':
M = 5
L = 0.35
n = 2
linking = "cluster-reverse"
inst_type = "cavemen"

i1, i2, jm = gen_special_jUFLP(n, M, L, linking, inst_type)

draw_jUFLP_inst(i1, i2, jm, "tmp/jUFLP_ex-new.dot")
save_inst(i1, i2, jm, "./instances/jUFLP-ex.json")

S, f, c, caves = i1

opt_order = UFLP_greedy_order(S)

N2s = N2RList(S)

D, nl = create_cover_DD(S, f, c, opt_order)

D.show(node_labels=nl)
print(f"order={opt_order}")
