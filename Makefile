# Aligning BDDs: the main makefile
# Contains all key instructions to reproduce the paper
# results (figures)
#
# (c) Alexey Bochkraev, Clemson University, 2020
# abochka@clemson.edu

######################################################################
## Files and directories
PREF=./
INST=$(PREF)/instances/raw
LOGS=$(PREF)/run_logs
ARC=./data_archive
FIGS=./figures
SOLVE=./solve_inst.py
BBLOG=./log_bb.py
PP=./post_processing
SCAL=./scal_test.py
STATS=./gen_lsizes_stats.py
CMP_LBs=./compare_simpl_LBs.py
PSCAL=./BB_orig_solve.py
######################################################################
## Numerical parameters
PAR_SOL=	4

### random dataset stats
LW_n=1000
LW_N=15
LW_Ps=0.3 0.6 0.9

### dataset generation parameters:
p=0.6# dataset generation parameter
n=10000# number of instances
N=15# number of variables per instance
n_LBs=1000

### scalability figure
SCAL_N=5 6 7 8 9 10 12 13 14 15 #16 17 18 19 20 #22 25 28
SCAL_K=50
SCAL_P=$(p)
SCAL_R=R
######################################################################
## calculated vars/parameters
SCAL_FILES=$(addsuffix .log, $(addprefix $(LOGS)/scal_$(SCAL_R),$(SCAL_N)))
LW_FILES_R=$(addsuffix .log, $(addprefix $(LOGS)/lwidths_Rp,$(LW_Ps)))
LW_FILES_N=$(addsuffix .log, $(addprefix $(LOGS)/lwidths_Np,$(LW_Ps)))

MAX_I=$(shell expr $(PAR_SOL) - 1 )
SOL_FILES_R = $(addsuffix .log, $(addprefix $(LOGS)/part.solved_R., $(shell seq -f %02g 0 $(MAX_I))))
SOL_FILES_N = $(addsuffix .log, $(addprefix $(LOGS)/part.solved_N., $(shell seq -f %02g 0 $(MAX_I))))

BB_FILES_R = $(addsuffix .log, $(addprefix $(LOGS)/part.BB_bounds_R., $(shell seq -f %02g 0 $(MAX_I))))
BB_FILES_N = $(addsuffix .log, $(addprefix $(LOGS)/part.BB_bounds_N., $(shell seq -f %02g 0 $(MAX_I))))

PSCAL_FILES = $(addsuffix .log, $(addprefix $(LOGS)/part.scal., $(shell seq -f %02g 0 $(MAX_I))))

DTE=$(shell date +%F)

MAX_I=$(shell expr $(PAR_SOL) - 1 )

.PHONY: all figures clean-raw-inst clean-insts move-logs figures/sample_BB_tree.png

######################################################################
## High-level recipes
.SECONDARY: # keep all the intermediary logfiles (will not work o/w)

all:

figures/sample_BB_tree.png: ./sample_BB_tree.py
	python ./sample_BB_tree.py -v -V 8 -n 10 -o ./run_logs/sample_BB_tree.dot && \
	dot -Tpng ./run_logs/sample_BB_tree.dot > ./figures/sample_BB_tree.png

figures: $(FIGS)/fig_sol_guessing_R.eps $(FIGS)/fig_sol_fireplace_R.eps $(FIGS)/fig_BB_gaps_R.eps $(FIGS)/LB.eps $(FIGS)/fig_sol_obj_hist_R.eps $(FIGS)/fig_sol_obj_int_R.eps $(FIGS)/guessing.tex

#figures_R figures_N

scalfig: $(FIGS)/fig_scal.eps

summary_figs: $(FIGS)/fig_summary_R.eps $(FIGS)/fig_summary_N.eps
	touch summary_figs

figures_R: $(FIGS)/fig_sol_guessing_R.eps $(FIGS)/fig_BB_gaps_R.eps $(FIGS)/fig_sol_fireplace_R.eps $(FIGS)/fig_sol_obj_hist_R.eps $(FIGS)/fig_sol_obj_int_R.eps
figures_N: $(FIGS)/fig_sol_guessing_N.eps $(FIGS)/fig_BB_gaps_N.eps $(FIGS)/fig_sol_fireplace_N.eps $(FIGS)/fig_sol_obj_hist_N.eps $(FIGS)/fig_sol_obj_int_N.eps

hists: $(FIGS)/fig_sol_obj_hist_R.eps $(FIGS)/fig_sol_obj_hist_N.eps
	touch hists

######################################################################
## Figure recipes
SOL_RECIPE = Rscript $(PP)/fig_$*.R -i $< -o $@
BB_RECIPE = Rscript $(PP)/fig_BB_$*.R -i $< -o $@

$(FIGS)/fig_sol_%_R.eps: $(LOGS)/solved_R.log $(PP)/fig_%.R
	$(SOL_RECIPE)

$(FIGS)/fig_sol_%_N.eps: $(LOGS)/solved_N.log $(PP)/fig_%.R
	$(SOL_RECIPE)

$(FIGS)/fig_BB_%_R.eps: $(LOGS)/BB_bounds_R.log $(PP)/fig_BB_%.R
	$(BB_RECIPE)

$(FIGS)/fig_BB_%_N.eps: $(LOGS)/BB_bounds_N.log $(PP)/fig_BB_%.R
	$(BB_RECIPE)

$(FIGS)/fig_scal.eps: $(LOGS)/scal_par.log $(PP)/fig_scal.R
	Rscript $(PP)/fig_scal.R -i $< -o $@

$(FIGS)/fig_summary_R.eps: DS_FLAG=R
$(FIGS)/fig_summary_N.eps: DS_FLAG=N

$(FIGS)/fig_summary_%.eps: $(LOGS)/lwidths_%.log $(PP)/fig_summary.R
	Rscript $(PP)/fig_summary.R -i $< -o $@

$(FIGS)/guessing.tex: $(LOGS)/solved_R.log $(PP)/tab_guessing.R
	Rscript $(PP)/tab_guessing.R -i $< -o $@

######################################################################
## scalability figure -- parallel implementation
$(LOGS)/scal_par.log: $(PSCAL_FILES)
	python $(PSCAL) --header > $@ && \
	tail -qn +2 $(PSCAL_FILES) >> $@ && \
	rm -rf $(PSCAL_FILES)

$(LOGS)/part.scal.%.log: $(INST)/scal/instances.list
	python $(PSCAL) -i $<.$* -o $@ -d $(INST)/scal/

$(INST)/scal/instances.list: ./gen_BDD_pair.py
	@echo Preparing instances for the scalability set...
	i=0 && \
	rm -rf $(INST)/scal && \
	mkdir -p $(INST)/scal && \
	for S in $(SCAL_N); do \
		python ./gen_BDD_pair.py -s $$(( ${SCAL_K} * $${i} )) -n $(SCAL_K) -v $$S -p $p -RU $(INST)/scal/ > /dev/null &&\
		i=$$(( $${i} + 1 ));\
	done && \
	ls $(INST)/scal | grep -Po "A\\K[^\\.]*" > $@ && \
	split -d -nl/$(PAR_SOL) $@ $(INST)/scal/instances.list.

# was: ls $(INST)/$*/A*.bdd | grep -Po "$(INST)/$*/A\\K[^\\.]*" > $@ &&\
# basename -a `ls $(INST)/$*/A*.bdd` | sed -n -e's/^A//p' | sed -n -e's/\.bdd//p' > $@ && \
#
######################################################################
## Main calculations (creating .log-files)
.SECONDEXPANSION:
$(LOGS)/solved_%.log: $$(SOL_FILES_%)
	python $(SOLVE) --header > $@ && \
	tail -qn +2 $(SOL_FILES_$*) >> $@


$(LOGS)/part.solved_R.%.log: $(INST)/R/instances.list $(SOLVE)
	python $(SOLVE) -i $<.$* -o $@ -d $(INST)/R/
######################################################################
$(FIGS)/LB.eps: $(LOGS)/LBs.log $(PP)/fig_LBs.R
	Rscript $(PP)/fig_LBs.R -i $< -o $@

$(LOGS)/LBs.log: $(CMP_LBs)
	mkdir -p $(INST)/LBs && \
	python ./gen_BDD_pair.py -n $(n_LBs) -v $N -p $p -RU $(INST)/LBs/ > $(LOGS)/$(DTE)_gen_$(N)var_LBs.log && \
	ls $(INST)/LBs | grep -Po "A\\K[^\\.]*" > $(INST)/LBs/instances.list && \
	python $(CMP_LBs) -d $(INST)/LBs/ -l $(INST)/LBs/instances.list > $@

## Generating and solving instances
$(INST)/%/instances.list:
	@echo Preparing the $*-instances dataset...
	mkdir -p $(INST)/$* && \
	if [ -f $(ARC)/dataset_$*.tar.gz ]; then tar -zxmf $(ARC)/dataset_$*.tar.gz; \
	else \
	python ./gen_BDD_pair.py -n $n -v $N -p $p -$*U $(INST)/$*/ > $(LOGS)/$(DTE)_gen_$(N)var_R.log; fi && \
	ls $(INST)/$* | grep -Po "A\\K[^\\.]*" > $@ && \
	split -d -nl/$(PAR_SOL) $@ $(INST)/$*/instances.list.

# was: ls $(INST)/$*/A*.bdd | grep -Po "$(INST)/$*/A\\K[^\\.]*" > $@ &&\
# basename -a `ls $(INST)/$*/A*.bdd` | sed -n -e's/^A//p' | sed -n -e's/\.bdd//p' > $@ && \
#
######################################################################
## Main calculations (creating .log-files)
.SECONDEXPANSION:
$(LOGS)/solved_%.log: $$(SOL_FILES_%)
	python $(SOLVE) --header > $@ && \
	tail -qn +2 $(SOL_FILES_$*) >> $@


$(LOGS)/part.solved_R.%.log: $(INST)/R/instances.list $(SOLVE)
	python $(SOLVE) -i $<.$* -o $@ -d $(INST)/R/

$(LOGS)/part.solved_N.%.log: $(INST)/N/instances.list $(SOLVE)
	python $(SOLVE) -i $<.$* -o $@ -d $(INST)/N/

$(LOGS)/BB_bounds_%.log: $$(BB_FILES_%)
	python $(BBLOG) --header > $@ && \
	tail -qn +2 $(BB_FILES_$*) >> $@

$(LOGS)/part.BB_bounds_R.%.log: $(INST)/R/instances.list $(SOLVE)
	python $(BBLOG) -i $<.$* -o $@ -d $(INST)/R/

$(LOGS)/part.BB_bounds_N.%.log: $(INST)/N/instances.list $(SOLVE)
	python $(BBLOG) -i $<.$* -o $@ -d $(INST)/N/

## scalability experiment
$(LOGS)/scal_$(SCAL_R)%.log: $(SCAL)
	mkdir -p $(INST)/sc$* && \
	python ./gen_BDD_pair.py -n $(SCAL_K) -v $* -p $(SCAL_P) -U$(SCAL_R) $(INST)/sc$* > $(INST)/sc$*/gen.log && \
	basename -a -s .bdd $(INST)/sc$*/A*.bdd | sed 's/A//' > $(INST)/sc$*/inst.list && \
	python $(SCAL) -l $(INST)/sc$*/inst.list -d $(INST)/sc$*/ > $@

$(LOGS)/scalability.log: $(SCAL) $(SCAL_FILES)
	python $(SCAL) --header > $@ && \
	cat $(SCAL_FILES) >> $@ && \
	tar --remove-files -czf $(ARC)/scalability_p$(SCAL_P)$(SCAL_R)x$(SCAL_K).tar.gz $(INST)/sc* && \
	rm -rf $(INST)/sc$* && \
	rm $(LOGS)/scal_N*.log
## end of scalability experiment

## random dataset stats
GEN_LW_LOGS =	mkdir -p $(INST)/ds_stats/$(DS_FLAG)p$* && \
	python ./gen_BDD_pair.py -n $(LW_n) -v $(LW_N) -p $* -$(DS_FLAG)U $(INST)/ds_stats/$(DS_FLAG)p$* > /dev/null && \
	ls $(INST)/ds_stats/$(DS_FLAG)p$*/*.bdd > $(INST)/ds_stats/$(DS_FLAG)p$*/lwidths_$(DS_FLAG)p$*.list && \
	python $(STATS) -s $* $(INST)/ds_stats/$(DS_FLAG)p$*/lwidths_$(DS_FLAG)p$*.list > $@ && \
	rm -rf $(INST)/ds_stats/$(DS_FLAG)p$*

$(LOGS)/lwidths_Np%.log: $(STATS)
	$(GEN_LW_LOGS)

$(LOGS)/lwidths_Rp%.log: $(STATS)
	$(GEN_LW_LOGS)

.SECONDEXPANSION:
$(LOGS)/lwidths_%.log: $$(LW_FILES_%)
	ls $(LOGS)/lwidths_$(DS_FLAG)p*.log && \
	tail -qn +2 $(LOGS)/lwidths_$(DS_FLAG)p*.log > $@

## end of random dataset stats

######################################################################
# auxiliary recipes

prep_dirs:
	mkdir -p $(INST)/R $(INST)/N
	mkdir -p $(LOGS)

check_src:
	egrep -nr --color 'TODO|FIXME|BUG|NOTE'

install_R_pkg:
	Rscript ./aux/R_install_packages.R

move-logs:
	@echo moving logs away from $(LOGS) to $(ARC)
	rm -f $(LOGS)/part.* && \
	tar --remove-files -czf $(ARC)/$(DTE)-logs.tar.gz -C $(LOGS) .

# clean recipes
clean-raw-inst:
	@echo Cleaning up raw instances in $(INST)...
	rm -f $(INST)/R/*.bdd
	rm -f $(INST)/R/*.list*
	rm -f $(INST)/N/*.bdd
	rm -f $(INST)/N/*.list*
	rm -rf $(INST)/ds_stats
	rm -f $(INST)/*.bdd

clean-archive:
	@echo Cleaning up instance archives...
	rm -f $(ARC)/*.tar.gz

clean-logs:
	@echo Cleaning up log files...
	rm -f $(LOGS)/*.log

clean-figures:
	@echo Cleaning up figures...
	rm -f $(FIGS)/*.eps

clean: clean-raw-inst clean-logs clean-figures clean-archive
clean-tmp: clean-raw-inst clean-logs

test: $(LOGS)/BB.test
	@echo Running tests...

$(LOGS)/BB.test: $(LOGS)/BB_bounds_R.log $(LOGS)/solved_R.log tests/BB_log_correct.R
	Rscript tests/BB_log_correct.R -b $(LOGS)/BB_bounds_R.log -s $(LOGS)/solved_R.log > $(LOGS)/BB.test
