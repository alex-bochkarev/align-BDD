# Aligning BDDs: the main makefile
# Contains all key instructions to reproduce the paper
# results (figures)
#
# (c) Alexey Bochkraev, Clemson University, 2020
# abochka@clemson.edu

# Makefile todo-list:
# TODO add random seed / other means for the dataset reproducibilty
# TODO add flexible calculation root directory
# TODO add parallel instances processing (if not yet)
# TODO make summary/stats figures for the dataset (# nodes per layer)
# TODO remove instances.list after the calculation

######################################################################
## Files and directories
INST=./instances/raw
LOGS=./run_logs
FIGS=./figures
OUT=./experiments
SOLVE=./solve_inst.py
BBLOG=./log_bb.py
PP=./post_processing
SCAL=./scal_test.py
STATS=./gen_lsizes_stats.py
STATS_LIST_R=./instances/raw/reduced/stats.list

######################################################################
## Numerical parameters

### scalability figure
SCAL_N=5 6 7 8 9 10 12 13 14 15 17 18
SCAL_K=5
SCAL_P=0.7

### random dataset stats
LW_n=50
LW_N=15
LW_Ps=0.3 0.5 0.7

### dataset generation parameters:
p=0.6# dataset generation parameter
n=100# number of instances
N=10# number of variables per instance

######################################################################
## calculated vars/parameters
SCAL_FILES=$(addsuffix .log, $(addprefix $(LOGS)/scal_N,$(SCAL_N)))
LW_FILES=$(addsuffix .log, $(addprefix $(LOGS)/lwidths_p,$(LW_Ps)))

DTE=$(shell date +%F)

.PHONY: all figures clean-raw-inst clean-insts

######################################################################
## High-level recipes

all: #instances/dataset_R.tar.gz
figures: $(FIGS)/fig_guessing_R.eps $(FIGS)/fig_BB_gaps_R.eps $(FIGS)/fig_fireplace_R.eps $(FIGS)/fig_obj_hist_R.eps $(FIGS)/fig_obj_int_R.eps $(FIGS)/fig_scal.eps

######################################################################
## Figure recipes

$(FIGS)/fig%_R.eps: $(LOGS)/solved_R.log $(PP)/fig%.R
	Rscript $(PP)/fig$*.R -i $< -o $@

$(FIGS)/fig_BB_%_R.eps: $(LOGS)/BB_bounds_R.log $(PP)/fig_BB_%.R
	Rscript $(PP)/fig_BB_$*.R -i $< -o $@

$(FIGS)/fig_scal.eps: $(LOGS)/scalability.log $(PP)/fig_scal.R
	Rscript $(PP)/fig_scal.R -i $< -o $@

$(FIGS)/fig_summary_R.eps: $(LOGS)/lwidths_R.log $(PP)/fig_summary.R
	Rscript $(PP)/fig_summary.R -i $< -o $@

## Generating and solving instances
$(INST)/reduced/instances.list:
	@echo Preparing the reduced instances dataset...
	test -f ./dataset_R.tar.gz && \
	tar -zxmf ./dataset_R.tar.gz \
	|| python ./gen_BDD_pair.py -n $n -v $N -p $p -RU $(INST)/reduced/ > $(LOGS)/$(DTE)_gen_$(N)var_R.log; \
	ls $(INST)/reduced/A*.bdd | grep -Po "$(INST)/reduced/A\\K[^\\.]*" > $(INST)/reduced/instances.list

instances/dataset_R.tar.gz: ./gen_BDD_pair.py $(INST)/reduced/*.bdd
	tar -czf ./instances/dataset_R.tar.gz $(INST)/*

######################################################################
## Main calculations (creating .log-files)

$(LOGS)/solved_R.log: $(INST)/reduced/instances.list $(SOLVE)
	python $(SOLVE) -i $< -o $@ -d $(INST)/reduced/

$(LOGS)/BB_bounds_R.log: $(INST)/reduced/instances.list $(BBLOG)
	python $(BBLOG) -i $< -o $@

## scalability experiment
$(LOGS)/scal_N%.log: $(SCAL)
	mkdir -p $(INST)/sc$* && \
	python ./gen_BDD_pair.py -n $(SCAL_K) -v $* -p $(SCAL_P) -U $(INST)/sc$* > $(INST)/sc$*/gen.log && \
	basename -a -s .bdd $(INST)/sc$*/A*.bdd | sed 's/A//' > $(INST)/sc$*/inst.list && \
	python $(SCAL) -l $(INST)/sc$*/inst.list -d $(INST)/sc$*/ > $@

$(LOGS)/scalability.log: $(SCAL) $(SCAL_FILES)
	python $(SCAL) --header > $@ && \
	cat $(wildcard $(LOGS)/scal_N*.log) >> $@ && \
	tar --remove-files -czf ./instances/scalability_p$(SCAL_P)x$(SCAL_K).tar.gz $(INST)/sc* && \
	rm -rf $(INST)/sc$* && \
	rm $(LOGS)/scal_N*.log
## end of scalability experiment

## random dataset stats
$(LOGS)/lwidths_p%.log: $(STATS)
	mkdir -p $(INST)/ds_stats/p$* && \
	python ./gen_BDD_pair.py -n $(LW_n) -v $(LW_N) -p $* -RU $(INST)/ds_stats/p$* > /dev/null && \
	ls $(INST)/ds_stats/p$*/*.bdd > $(INST)/ds_stats/p$*/lwidths_p$*.list && \
	python $(STATS) -s $* $(INST)/ds_stats/p$*/lwidths_p$*.list > $@ && \
	rm -rf $(INST)/ds_stats/p$*

$(LOGS)/lwidths_R.log: $(LW_FILES)
	tail -qn +2 $(LOGS)/lwidths_p*.log > $@ && \
	rm -rf $(INST)/ds_stats
## end of random dataset stats

######################################################################
# auxiliary recipes
check_src:
	egrep -nr --color 'TODO|FIXME|BUG|NOTE'

# clean recipes
clean-raw-inst:
	@echo Cleaning up raw instances in $(INST)...
	rm -f $(INST)/reduced/*.bdd
	rm -f $(INST)/reduced/*.list
	rm -f $(INST)/nonreduced/*.bdd
	rm -f $(INST)/nonreduced/*.list
	rm -f $(INST)/*.bdd

clean-inst:
	@echo Cleaning up instance archives...
	rm -f $(INST)/*.tar.gz

clean-logs:
	@echo Cleaning up log files...
	rm -f $(LOGS)/*.log

clean-figures:
	@echo Cleaning up figures...
	rm -f $(FIGS)/*.eps

clean: clean-raw-inst clean-inst clean-logs clean-figures
