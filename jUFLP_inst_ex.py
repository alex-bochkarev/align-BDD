"""A simple illustration for the special instance type for joint UFLP."""
from UFLP_2_cav import gen_special_jUFLP, draw_jUFLP_inst


# if __name__ == '__main__':
M = 5
L = 0.35
n = 2
linking = "cluster-reverse"
inst_type = "cavemen"

i1, i2, jm = gen_special_jUFLP(n, M, L, linking, inst_type)

draw_jUFLP_inst(i1, i2, jm)
