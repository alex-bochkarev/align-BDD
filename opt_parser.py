"""An auxiliary script: parser for working with datasets of saved align-varseq
instances.
"""

######################################################################
## A parser for the generated optimal alignments
import re
from varseq import *
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

class alignments:
    def __init__(self, filename):
        self.fname = filename;
        ## fields to be read from the file
        self.instance_no = []
        self.perm_no = []
        self.As = []
        self.Bs = []
        self.As_aligned = []
        self.Bs_aligned = []
        self.no_opts = []

    def fmt_err(self, msg, ln):
        print("File format error: <{}, line {}> -- {}".format(self.fname, ln, msg))

    def get_instance(self,inst):
        if inst<0:
            inst = self.instance_no[inst]
        return [self.As[inst], self.Bs[inst]]

    def show_instance(self,inst):
        if inst<0:
            inst = self.instance_no[inst]

        print("| Ins/opt  | A                       | B                       | A-layers                | B-layers                |")
        print("+----------+-------------------------+-------------------------+-------------------------+-------------------------+")
        print('| {:<9}| {:<24}| {:<24}| {:<24}| {:<24}|'.format(self.instance_no[inst], str(self.As[inst].layer_var), str(self.Bs[inst].layer_var),str(self.As[inst].n), str(self.Bs[inst].n)))
        print("+----------+-------------------------+-------------------------+-------------------------+-------------------------+")
        for i in range(self.no_opts[inst]):
            print('| {:<9}| {:<24}| {:<24}| {:<24}| {:<24}|'.format("opt. #{}".format(i), str(self.As_aligned[inst][i].layer_var), str(self.Bs_aligned[inst][i].layer_var),str(self.As_aligned[inst][i].n), str(self.Bs_aligned[inst][i].n)))

        print("+----------+-------------------------+-------------------------+-------------------------+-------------------------+")

        print("Initial size: {}, optimal aligned size: {}".format(self.As[inst].size()+self.Bs[inst].size(), self.As_aligned[inst][0].size()+self.Bs_aligned[inst][0].size()))

    def fmt_instance(self,A,B):
        optsno, Aas,Bas = A.OA_bruteforce(B)

        print("| Ins/opt  | A                       | B                       | A-layers                | B-layers                |")
        print("+----------+-------------------------+-------------------------+-------------------------+-------------------------+")
        print('| {:<9}| {:<24}| {:<24}| {:<24}| {:<24}|'.format("adhoc", str(A.layer_var), str(B.layer_var),str(A.n), str(B.n)))
        print("+----------+-------------------------+-------------------------+-------------------------+-------------------------+")
        for i in range(optsno):
            print('| {:<9}| {:<24}| {:<24}| {:<24}| {:<24}|'.format("opt. #{}".format(i), str(Aas[i].layer_var), str(Bas[i].layer_var),str(Aas[i].n), str(Bas[i].n)))

        print("+----------+-------------------------+-------------------------+-------------------------+-------------------------+")

        print("Initial size: {}, optimal aligned size: {}".format(A.size()+B.size(), Aas[0].size()+Bas[0].size()))
    ## loads a log-file with generated optimal
    ## alignments
    def load(self, verbose = True):
        # prepare some regexps
        p_inst = re.compile("# instance no: ([0-9]+)")
        p_perm = re.compile("## Perm. no: ([0-9]+)\n")
        p_vars = re.compile(r'Vars: \[([0-9\s]+)\]\n')
        p_n = re.compile(r'n\s+: \[([0-9\s]+)\] \(sz=([0-9]+)\)\n')
        p_opts = re.compile(r'Opt alignments: ([0-9]+)\n')
        p_totsize = re.compile(r'\*\*Total size:\*\* ([0-9]+)\n')

        ln = 0
        with open(self.fname,"r") as f:
            for line in f:
                ln += 1
                m = p_inst.match(line)
                if m is None: continue # if it is not an instance description start -- skip

                self.instance_no.append(int(m[1]))
                m = p_perm.match(f.readline())
                ln+=1
                if m is None:
                    self.fmt_err("No perm no. entry after instance entry", ln)
                    return False
                self.perm_no.append(int(m[1]))

                if f.readline() == "A:\n":
                    mv = p_vars.match(f.readline())
                    mn = p_n.match(f.readline())
                    self.As.append(VarSeq(layer_vars = list(map(lambda x: int(x), mv[1].split())),
                                     layer_sizes = list(map(lambda x: int(x), mn[1].split()))))
                else:
                    self.fmt_err("Entry for A: expected", ln)
                    return False

                if f.readline() == "B:\n":
                    mv = p_vars.match(f.readline())
                    mn = p_n.match(f.readline())
                    self.Bs.append(VarSeq(layer_vars = list(map(lambda x: int(x), mv[1].split())),
                                     layer_sizes = list(map(lambda x: int(x), mn[1].split()))))
                else:
                    self.fmt_err("Entry for B: expected", ln)
                    return False

                ln += 7

                ## cross-check
                m = p_totsize.match(f.readline())

                if m is None or int(m[1]) != self.As[-1].size()+self.Bs[-1].size():
                    self.fmt_err("Total size is wrong ({} read, but {} expected)".format(int(m[1]), self.As[-1].size()+self.Bs[-1].size()), ln)
                    return False

                ln+=1
                m = p_opts.match(f.readline())

                if m is None:
                    self.fmt_err("No. of opt alignments entry is expected", ln)
                    return False

                self.no_opts.append(int(m[1]))
                ln+=1

                self.As_aligned.append([])
                self.Bs_aligned.append([])

                for al in range(self.no_opts[-1]):
                    _ = f.readline()
                    _ = f.readline()
                    _ = f.readline()

                    ## read A
                    mv = p_vars.match(f.readline())
                    mn = p_n.match(f.readline())

                    if mv is None or mn is None:
                        self.fmt_err("Error reading A-aligned, opt #{}".format(al),ln)
                        return False

                    self.As_aligned[-1].append(VarSeq(layer_vars = list(map(lambda x: int(x), mv[1].split())),
                                                 layer_sizes = list(map(lambda x: int(x), mn[1].split()))))
                    _ = f.readline()
                    ## read B
                    mv = p_vars.match(f.readline())
                    mn = p_n.match(f.readline())
                    if mv is None or mn is None:
                        self.fmt_err("Error reading B-aligned, opt #{}".format(al),ln)
                        return False

                    self.Bs_aligned[-1].append(VarSeq(layer_vars = list(map(lambda x: int(x), mv[1].split())),
                                                 layer_sizes = list(map(lambda x: int(x), mn[1].split()))))
                if len(self.instance_no) % 100 == 0 and verbose:
                    print("=", end="",flush=True)

        if verbose: print("\n{} instances loaded".format(len(self.instance_no)))
        return True

    # checks conditions / hypotheses over loaded instances
    def check_conditions(self, conditions, names, instances=-1):
        df = pd.DataFrame(data = {'Proposition':names,'∀ opt: holds':[False for i in range(len(names))], '∃ opt: holds':[True for i in range(len(names))], '% opt: holds':[0.0 for i in range(len(names))], 'Fails':[[] for i in range(len(names))],'No-opt-fails':[[] for i in range(len(names))]})

        opt_share = np.zeros(len(conditions))

        opt_all = [True for i in range(len(conditions))]
        opts = 0

        if instances == -1:
            I = len(self.instance_no)
        else:
            I = instances

        for i in range(I):
            opt_ex = [False for i in range(len(conditions))]
            for opt in range(self.no_opts[i]):
                opts += 1
                for c in range(len(conditions)):
                    if conditions[c](self.As[i],self.Bs[i],self.As_aligned[i][opt], self.Bs_aligned[i][opt]):
                        opt_ex[c] = True
                        opt_share[c]+=1
                    else:
                        opt_all[c] = False
                        df['Fails'][c].append([i, opt])

                if opts % 100 == 0:
                    print("=",end="", flush=True)
            df['∃ opt: holds'] = df['∃ opt: holds'] & opt_ex

            # save instances where a condition failed for *all* optimal alignments
            for c in range(len(conditions)):
                if not opt_ex[c]:
                    df['No-opt-fails'][c].append(i)

        df['∀ opt: holds'] = opt_all
        df['% opt: holds'] = opt_share / opts

        print("done")
        return df

    ## API to save a case
    def save_case(self, A,B,i_no = -1):
        with open(self.fname,"a") as outf:

            if i_no == -1:
                i = len(self.As)
            else:
                i = i_no

            outf.write("\n# instance no: {}\n".format(i))
            outf.write("## Perm. no: 0\n")
            res = A.OA_bruteforce(B)
            outf.write("A:\n{}\nB:\n{}\n**Total size:** {}\nOpt alignments: {}\n".format(A,B,A.size()+B.size(), res[0]))
            outf.write("----------------------------------------------\n")
            for o in range(res[0]):
                outf.write("### Optimal alignment no: {}\n".format(o))
                outf.write("**A(aligned):**\n{}\n**B(aligned):**\n{}\n".format(res[1][o],res[2][o]))
                outf.write("----------------------------------------------\n")
            outf.write("==============================================\n")
        # save the instance into the currently loaded data-structure
        self.As.append(A)
        self.Bs.append(B)
        self.As_aligned.append(res[1])
        self.Bs_aligned.append(res[2])
        self.no_opts.append(res[0])
        self.instance_no.append(i)
        self.perm_no.append(-1)
# ## testing code

# parser = alignments("./opt_als.md")
# parser.load()

# sns.distplot(parser.no_opts,rug=True, kde=False)
# plt.xlabel("No. of optimal alignments")
# plt.ylabel("No. of 7-var instances (out of {})".format(len(parser.instance_no)))
# plt.show()

# for i in [5,14,23]:
#     parser.show_instance(i)
