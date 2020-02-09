"""
Benchmarks several heuristic algos against a given set
of varseq instances (from a given folder)

(c) A. Bochkarev, Clemson University, 2019
abochka@clemson.edu
"""
import sys
import requests # needed for Telegram
import numpy as np

import opt_parser as als
import varseq as vs
from glob import glob
from copy import deepcopy, copy

# for debugging // profiling
import cProfile
import pstats


def tg_sendtext(bot_message):
    bot_token = '<token>'
    bot_chatID = '<chatid>'
    send_text = 'https://api.telegram.org/bot' + bot_token + '/sendMessage?chat_id=' + bot_chatID + '&parse_mode=Markdown&text=' + bot_message
    response = requests.get(send_text)
    return response.json()

def toA(A,B):
    Bp = B.align_to(A)
    return [A.size()+Bp.size(),A.layer_var]

def toB(A,B):
    return toA(B,A)

# min { toA, toB}
def minAB(A,B):
    sA,vA = toA(A,B)
    sB,vB = toA(B,A)

    if sA<=sB:
        return sA,vA
    else:
        return sB,vB

def simpl_greedy_swaps(A,B):
    """finds a good alignment target with greedy swaps

    Makes swaps while it can improve the objective
    (picking the best improvement at each step)
    """
    cur_size, T = minAB(A,B)
    T = vs.VarSeq(T,[1 for i in range(len(T))])
    improvement_possible = True

    while improvement_possible:
        improvement_possible = False
        Tp = deepcopy(T)

        for i in range(len(T)-1):
            a = T.layer_var[i];
            Ap = A.align_to(T.slide(a,i-1))
            Bp = B.align_to(T.slide(a,i-1))

            if Ap.size() + Bp.size() < cur_size:
                cur_size = Ap.size() + Bp.size()
                Tp = T.slide(a,i-1)
                improvement_possible = True
        T = Tp

    return [A.align_to(T).size() + B.align_to(T).size(),T]

# gr sifts
def simpl_greedy_sifts(A,B,passes=-1):
    """finds a good alignment target with greedy sifts

    Makes sifts while it can improve the objective
    or for a given number of passes,
    picking the best improvement at each step.
    (one pass = examining all the variables)
    """
    cur_size, T = minAB(A,B)
    T = vs.VarSeq(T,[1 for i in range(len(T))])
    improvement_possible = True
    pass_no = 0

    while improvement_possible:
        improvement_possible = False
        pass_no += 1
        Tp = deepcopy(T)

        for i in range(len(T)):
            a = T.layer_var[i]
            for j in range(len(T)):
                Ap = A.align_to(T.slide(a,j))
                Bp = B.align_to(T.slide(a,j))

                if Ap.size() + Bp.size() < cur_size:
                    cur_size = Ap.size() + Bp.size()
                    Tp = T.slide(a,j)
                    improvement_possible = True
        T = Tp
        if pass_no == passes:
            break

    return [A.align_to(T).size() + B.align_to(T).size(),T]

def simpl_gsifts_1p(A,B):
    return simpl_greedy_sifts(A,B,passes=1)

def simpl_gsifts_2p(A,B):
    return simpl_greedy_sifts(A,B,passes=2)

def simpl_gsifts_3p(A,B):
    return simpl_greedy_sifts(A,B,passes=3)

def simpl_5random(A,B):
    best_size = -1
    for i in range(5):
        o = np.random.permutation(A.layer_var)
        Ap = A.align_to(o)
        Bp = B.align_to(o)
        if best_size<0 or Ap.size()+Bp.size()<best_size:
            best_size = Ap.size()+Bp.size()
            best_o = copy(Ap.layer_var)
    return [best_size, best_o]

def simpl_greedy_2sifts(A,B, passes=-1):
    """finds a good alignment target with greedy pairs of sifts

    Makes sifts while it can improve the objective,
    or for a given number of passes,
    picking the best improvement at each step.
    (one pass = examining all the variables.)
    """
    cur_size, T = minAB(A,B)
    T = vs.VarSeq(T,[1 for i in range(len(T))])
    improvement_possible = True
    pass_no = 0

    while improvement_possible:
        improvement_possible = False
        pass_no += 1
        Tp = deepcopy(T)

        for i1 in range(len(T)):
            a1 = T.layer_var[i1]
            for j1 in range(len(T)):
                for i2 in range(len(T)):
                    a2 = T.layer_var[i2]
                    for j2 in range(len(T)):
                        Ap = A.align_to(T.slide(a1,j1).slide(a2,j2))
                        Bp = B.align_to(T.slide(a1,j1).slide(a2,j2))

                        if Ap.size() + Bp.size() < cur_size:
                            cur_size = Ap.size() + Bp.size()
                            Tp = T.slide(a1,j1).slide(a2,j2)
                            improvement_possible = True
        T = Tp
        if pass_no == passes:
            break

    return [A.align_to(T).size() + B.align_to(T).size(),T]

def simpl_g2sifts_1p(A,B):
    return simpl_greedy_2sifts(A,B,1)

def simpl_g2sifts_2p(A,B):
    return simpl_greedy_2sifts(A,B,2)

######################################################################
## an optimized version
# TODO: invent
def fastslide(X,a,p):
    if X[p]==a:
        return X

    X = [x for x in X]
    # so we definitely change something

    for i in range(len(X)):
        if X[i] != a and p > i:
            continue # nothing to do, yet

        elif X[i] == a:
            # so, p>i -- slide down
            j = i
            while j < p:
                X[j] = X[j+1]
                j += 1

            X[p] = a
            break

        elif p == i:
            # so, p(a)>i -- slide up
            buf1 = X[i]
            X[i] = a
            j = i+1
            while X[j]!=a:
                buf2 = X[j]
                X[j] = buf1
                buf1 = buf2
                j += 1

            X[j] = buf1
            break

    return X

def fast_greedy_sifts(A,B):
    """finds a good alignment target with greedy sifts

    Makes sifts while it can improve the objective
    (picking the best improvement at each step)
    """
    cur_size, T = minAB(A,B)
    improvement_possible = True
    no_of_ops = 0

    while improvement_possible:
        improvement_possible = False
        Tp = T

        for i in range(len(T)):
            a = T[i]
            for j in range(len(T)):
                Tpp = fastslide(T,a,j)
                Ap = A.align_to(Tpp)
                Bp = B.align_to(Tpp)

                no_of_ops += 2*abs(i-j) + 2*len(A)

                cand_size = Ap.size() + Bp.size()
                if cand_size < cur_size:
                    cur_size = cand_size
                    Tp = Tpp
                    no_of_ops += abs(i-j)
                    improvement_possible = True
        T = Tp

    return [T, A.align_to(T).size() + B.align_to(T).size(),no_of_ops]

def fast_greedy_2sifts(A,B):
    """finds a good alignment target with greedy pairs of sifts

    Makes sifts while it can improve the objective
    (picking the best improvement at each step)
    """
    cur_size, T  = minAB(A,B)
    improvement_possible = True
    no_of_ops = 0

    while improvement_possible:
        # print("obj:{}".format(cur_size))
        improvement_possible = False
        Tp = T
        for i1 in range(len(T)-1):
            a1 = T[i1]
            for j1 in range(i1+1,len(T)):
                Tpi = fastslide(T,a1,j1)

                for i2 in range(i1+1,len(T)-1):
                    a2 = T[i2]
                    for j2 in range(i2+1,len(T)):
                        Tpp = fastslide(Tpi,a2,j2)
                        Ap = A.align_to(Tpp)
                        Bp = B.align_to(Tpp)
                        no_of_ops += 2*len(A)
                        cand_size = Ap.size() + Bp.size()
                        if cand_size < cur_size:
                            cur_size = cand_size
                            Tp = Tpp
                            no_of_ops += abs(i1-j1)+abs(i2-j2)
                            improvement_possible = True
        T = Tp

    return [T, A.align_to(T).size() + B.align_to(T).size(),no_of_ops]


# heuristics = [
#     minAB,
#     greedy_swaps,
#     greedy_sifts,
#     greedy_2sifts
# ]

# heu_description = [
#     "minAB",
#     "gswaps",
#     "gsifts",
#     "g2sifts",
#     "smart_last_elem"
# ]

def orig_simpl(A,B,simpl):
    Ap = deepcopy(A); Bp = deepcopy(B)
    Ap.align_to(simpl["simpl_BB"][2])
    Bp.align_to(simpl["simpl_BB"][2])
    return [Ap.size()+Bp.size(),simpl["simpl_BB"][1], simpl["simpl_BB"][2]]

def orig_gsifts1p(A,B,simpl):
    ApA = deepcopy(A);BpA = deepcopy(B)
    ApB = deepcopy(A);BpB = deepcopy(B)
    ApB.gsifts(BpB)
    BpA.gsifts(ApA)
    sA = ApA.size()+BpA.size()
    sB = ApB.size()+BpB.size()
    if sA < sB:
        return [ApA.size()+BpA.size(), 0, ApA.vars]
    else:
        return [ApB.size()+BpB.size(), 0, ApB.vars]

def orig_bestAB(A,B,simpl):
    Ap = deepcopy(A); Bp = deepcopy(B)
    Ap.align_to(B.vars,inplace=True)
    Bp.align_to(A.vars,inplace=True)
    if Ap.size()+B.size() <= A.size()+Bp.size():
        return [Ap.size()+B.size(),0,Ap.vars]
    else:
        return [A.size()+Bp.size(),0,Bp.vars]

def orig_5random(A,B,simpl):
    best_size = -1
    for i in range(5):
        Ap = deepcopy(A);
        Bp = deepcopy(B);
        o = np.random.permutation(A.vars)
        Ap.align_to(o,inplace=True)
        Bp.align_to(o,inplace=True)
        if best_size<0 or Ap.size()+Bp.size()<best_size:
            best_size = Ap.size()+Bp.size()
            best_o = copy(Ap.vars)
    return [best_size, 0, best_o]

# heuristic structure
# CODE - FUNC - LEGEND
SIMPL_HEU = [
    # ["simpl_toA",toA, "align to A"],
    # ["simpl_toB",toB, "align to B"],
    ["simpl_minAB",minAB, "best of A and B"],
    ["simpl_gswaps",simpl_greedy_swaps,"greedy swaps"],
    ["simpl_5random",simpl_5random,"best of 5 random orders"],
    ["simpl_gsifts_1p",simpl_gsifts_1p,"greedy sifts (1 pass)"],
    ["simpl_gsifts_2p",simpl_gsifts_2p,"greedy sifts (2 passes)"],
    ["simpl_gsifts_3p",simpl_gsifts_3p,"greedy sifts (3 passes)"]
    # ["simpl_g2sifts_1p",simpl_g2sifts_1p,"greedy 2-sifts (1 pass)"]
]

ORIG_HEU = [
    ["orig_simpl",orig_simpl,"simplified problem"],
    ["orig_gsifts1p",orig_gsifts1p,"greedy BDD sifts (1 pass)"],
    ["orig_bestAB",orig_bestAB,"best of A and B"],
    ["orig_5random",orig_5random,"best of 5 random orders"]
]
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("USAGE: {} <inst_directory>/".format(sys.argv[0]))
        exit()

    inst_dir = sys.argv[1]

    print("inst_no,opt_obj,"+"_obj,".join(heu_description)+"_obj"+ ","+"_steps,".join(heu_description)+"_steps")

    pr = cProfile.Profile()
    tg_sendtext("LOCAL JOB: heuristics testing started, instance dir: {}".format(inst_dir))
    n = 0
    pr.enable()

    for file in glob(inst_dir+"vs*.md"):
        inst_no = file.split('.')[-2]
        inst_no = inst_no.split('/vs')[-1]

        inst = als.alignments(file)
        inst.load(verbose=False)

        A = inst.As[0]
        B = inst.Bs[0]

        objs = []
        steps = []

        for h in heuristics:
            T, obj,st = h(A,B)
            objs.append(str(obj))
            steps.append(str(st))

        print(inst_no+","+str(inst.As_aligned[0][0].size()+inst.Bs_aligned[0][0].size())+","+",".join(objs)+","+",".join(steps))

        sys.stdout.flush()
        n+= 1

    pr.disable()
    ps = pstats.Stats(pr).sort_stats('time')

    ps.dump_stats("./run_profile.dmp")
    tg_sendtext("LOCAL (heuristics): done. {} files processed".format(n))

