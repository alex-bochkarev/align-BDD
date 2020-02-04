"""
Benchmarks several heuristic algos against a given set
of varseq instances (from a given folder)

(c) A. Bochkarev, Clemson University, 2019
abochka@clemson.edu
"""
import sys
import requests # needed for Telegram

import opt_parser as als
import varseq as vs
from glob import glob
from copy import deepcopy

# for debugging // profiling
import cProfile
import pstats

def tg_sendtext(bot_message):
    bot_token = '<token>'
    bot_chatID = '<chatid>'
    send_text = 'https://api.telegram.org/bot' + bot_token + '/sendMessage?chat_id=' + bot_chatID + '&parse_mode=Markdown&text=' + bot_message
    response = requests.get(send_text)
    return response.json()

# IDEA: implement a more careful steps counter within classes?

# min { toA, toB}
def minAB(A,B):
    Ap = A.align_to(B)
    Bp = B.align_to(A)

    toA = A.size() + Bp.size()
    toB = Ap.size() + B.size()
    if  toA < toB:
        return [A.layer_var, toA, 2*len(A)]
    else:
        return [B.layer_var, toB, 2*len(B)]

# gr swaps
def greedy_swaps(A,B):
    """finds a good alignment target with greedy swaps

    Makes swaps while it can improve the objective
    (picking the best improvement at each step)
    """
    T, cur_size, _ = minAB(A,B)
    T = vs.VarSeq(T,[1 for i in range(len(T))])
    improvement_possible = True

    no_of_ops = 0

    while improvement_possible:
        improvement_possible = False
        Tp = deepcopy(T)
        no_of_ops += len(Tp)

        for i in range(len(T)-1):
            a = T.layer_var[i];
            Ap = A.align_to(T.slide(a,i-1))
            Bp = B.align_to(T.slide(a,i-1))

            no_of_ops += 2*len(A)

            if Ap.size() + Bp.size() < cur_size:
                cur_size = Ap.size() + Bp.size()
                Tp = T.slide(a,i-1)
                no_of_ops += 1
                improvement_possible = True
        T = Tp

    return [T, A.align_to(T).size() + B.align_to(T).size(), no_of_ops]

# gr sifts
def greedy_sifts(A,B):
    """finds a good alignment target with greedy sifts

    Makes sifts while it can improve the objective
    (picking the best improvement at each step)
    """
    T, cur_size, _= minAB(A,B)
    T = vs.VarSeq(T,[1 for i in range(len(T))])
    improvement_possible = True
    no_of_ops = 0

    while improvement_possible:
        improvement_possible = False
        Tp = deepcopy(T)
        no_of_ops += len(Tp)

        for i in range(len(T)):
            a = T.layer_var[i]
            for j in range(len(T)):
                Ap = A.align_to(T.slide(a,j))
                Bp = B.align_to(T.slide(a,j))

                no_of_ops += 2*abs(i-j) + 2*len(A)

                if Ap.size() + Bp.size() < cur_size:
                    cur_size = Ap.size() + Bp.size()
                    Tp = T.slide(a,j)
                    no_of_ops += abs(i-j)
                    improvement_possible = True
        T = Tp

    return [T, A.align_to(T).size() + B.align_to(T).size(),no_of_ops]

# gr 2sifts
# TODO: speed up / way too slow to use
def greedy_2sifts(A,B):
    """finds a good alignment target with greedy pairs of sifts

    Makes sifts while it can improve the objective
    (picking the best improvement at each step)
    """
    T, cur_size, _ = minAB(A,B)
    T = vs.VarSeq(T,[1 for i in range(len(T))])
    improvement_possible = True
    no_of_ops = 0

    while improvement_possible:
        improvement_possible = False
        Tp = deepcopy(T)
        no_of_ops += len(Tp)

        for i1 in range(len(T)):
            a1 = T.layer_var[i1]
            for j1 in range(len(T)):
                for i2 in range(len(T)):
                    a2 = T.layer_var[i2]
                    for j2 in range(len(T)):
                        Ap = A.align_to(T.slide(a1,j1).slide(a2,j2))
                        Bp = B.align_to(T.slide(a1,j1).slide(a2,j2))
                        no_of_ops += 2*len(A)

                        if Ap.size() + Bp.size() < cur_size:
                            cur_size = Ap.size() + Bp.size()
                            Tp = T.slide(a1,j1).slide(a2,j2)
                            no_of_ops += abs(i1-j1)+abs(i2-j2)
                            improvement_possible = True
        T = Tp

    return [T, A.align_to(T).size() + B.align_to(T).size(),no_of_ops]

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
    T, cur_size, _= minAB(A,B)
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
    T, cur_size, _ = minAB(A,B)
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

######################################################################
# TODO: smart last-el greedy
def smart_last_elem(A,B):
    cur_last_pos = len(A)-1

    Ap = deepcopy(A)
    Bp = deepcopy(B)
    A_set = set()
    B_set = set()

    no_of_ops = 0

    for n in range(len(A)-1,0,-1):
        # n is the current ``last'' position to choose a var for
        cur_cost = -1
        i_cand = -1
        max_Bpos_seen = -1

        for i in range(n,-1,-1):
            if Bp.p[Ap.layer_var[i]] < max_Bpos_seen or Ap.layer_var[i] in A_set or Ap.layer_var[i] in B_set:
                no_of_ops += 1
                continue # this var is ``dominated'' -- can't be the last one

            cand_cost = Ap.S(Ap.layer_var[i],n)+Bp.S(Ap.layer_var[i],n)
            no_of_ops += abs(i-n)*2
            if cand_cost < cur_cost or cur_cost==-1:
                cur_cost = cand_cost
                i_cand = i

            if Bp.p[Ap.layer_var[i]] > max_Bpos_seen:
                max_Bpos_seen = Bp.p[Ap.layer_var[i]]

        A_covered = set(Ap.layer_var[i_cand:n])
        B_covered = set(Bp.layer_var[Bp.p[Ap.layer_var[i_cand]]:n])

        j_A = i_cand
        j_B = Bp.p[ Ap.layer_var[i_cand] ]

        Ap = Ap.slide(Ap.layer_var[i_cand],n); Bp = Bp.slide(Ap.layer_var[i_cand],n)
        no_of_ops += abs(i_cand - n) + abs(Bp.p[Ap.layer_var[i_cand]] - n)

        # now, reshuffle the covered elements
        # O(N) procedure

        if len(A_covered)>1:
            for i in range(n):
                no_of_ops += 1
                if Bp.layer_var[i] in A_covered:
                    Ap.layer_var[j_A] = Bp.layer_var[i]
                    Ap.p[Bp.layer_var[i]] = j_A
                    j_A += 1

                if Ap.layer_var[i] in B_covered:
                    Bp.layer_var[j_B] = Ap.layer_var[i]
                    Bp.p [ Ap.layer_var[i] ] = j_B
                    j_B += 1

        A_set.update(A_covered)
        B_set.update(A_covered)

    return [Ap.layer_var, Ap.size()+Bp.size(),no_of_ops]

heuristics = [
    minAB,
    greedy_swaps,
    greedy_sifts,
    greedy_2sifts,
    smart_last_elem
]

heu_description = [
    "minAB",
    "gswaps",
    "gsifts",
    "g2sifts",
    "smart_last_elem"
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
