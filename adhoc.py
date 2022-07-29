import json
from jUFLP_utils import load_inst, draw_jUFLP_inst

i1, i2, link = load_inst("./instances/jUFLP_cm/inst_1.json")

draw_jUFLP_inst(i1, i2, link, filename="tmp/jUFLP-check.dot")
