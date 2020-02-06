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
PAR_SOL=8

### random dataset stats
LW_n=50
LW_N=15
LW_Ps=0.3 0.5 0.7

### dataset generation parameters:
p=0.6# dataset generation parameter
n=100# number of instances
N=10# number of variables per instance

### scalability figure
SCAL_N=5 6 7 8 9 10 12 13 14 15 16 17 18 19 20
SCAL_K=5
SCAL_P=$(p)
SCAL_R=N
######################################################################
## calculated vars/parameters
SCAL_FILES=$(addsuffix .log, $(addprefix $(LOGS)/scal_$(SCAL_R),$(SCAL_N)))
LW_FILES_R=$(addsuffix .log, $(addprefix $(LOGS)/lwidths_Rp,$(LW_Ps)))
LW_FILES_N=$(addsuffix .log, $(addprefix $(LOGS)/lwidths_Np,$(LW_Ps)))

MAX_I=$(shell expr $(PAR_SOL) - 1 )
SOL_FILES_R = $(addsuffix .log, $(addprefix $(LOGS)/part.solved_R., $(shell seq -f %02g 0 $(MAX_I))))
SOL_FILES_N = $(addsuffix .log, $(addprefix $(LOGS)/part.solved_N., $(shell seq -f %02g 0 $(MAX_I))))

DTE=$(shell date +%F)

.PHONY: all figures clean-raw-inst clean-insts

######################################################################
## High-level recipes
.SECONDARY: # keep all the intermediary logfiles (will not work o/w)

all: #instances/dataset_R.tar.gz
figures: $(FIGS)/fig_sol_guessing_N.eps $(FIGS)/fig_BB_gaps_N.eps $(FIGS)/fig_sol_fireplace_N.eps $(FIGS)/fig_sol_obj_hist_N.eps $(FIGS)/fig_sol_obj_int_N.eps $(FIGS)/fig_scal.eps $(FIGS)/fig_sol_obj_hist_R.eps

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

$(FIGS)/fig_scal.eps: $(LOGS)/scalability.log $(PP)/fig_scal.R
	Rscript $(PP)/fig_scal.R -i $< -o $@

$(FIGS)/fig_summary_R.eps: DS_FLAG=R
$(FIGS)/fig_summary_N.eps: DS_FLAG=N

$(FIGS)/fig_summary_%.eps: $(LOGS)/lwidths_%.log $(PP)/fig_summary.R
	Rscript $(PP)/fig_summary.R -i $< -o $@

## Generating and solving instances
$(INST)/%/instances.list:
	@echo Preparing the $*-instances dataset...
	mkdir -p $(INST)/$* && \
	if [ -f $(INST)/dataset_$*.tar.gz ]; then tar -zxmf $(INST)/dataset_$*.tar.gz; \
	else \
	python ./gen_BDD_pair.py -n $n -v $N -p $p -$*U $(INST)/$*/ > $(LOGS)/$(DTE)_gen_$(N)var_R.log; fi && \
	ls $(INST)/$*/A*.bdd | grep -Po "$(INST)/$*/A\\K[^\\.]*" > $@ &&\
	split -d -nl/$(PAR_SOL) $@ $(INST)/$*/instances.list.

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

$(LOGS)/BB_bounds_%.log: $(INST)/%/instances.list $(BBLOG)
	python $(BBLOG) -d $(INST)/$*/ -i $< -o $@

## scalability experiment
$(LOGS)/scal_$(SCAL_R)%.log: $(SCAL)
	mkdir -p $(INST)/sc$* && \
	python ./gen_BDD_pair.py -n $(SCAL_K) -v $* -p $(SCAL_P) -U$(SCAL_R) $(INST)/sc$* > $(INST)/sc$*/gen.log && \
	basename -a -s .bdd $(INST)/sc$*/A*.bdd | sed 's/A//' > $(INST)/sc$*/inst.list && \
	python $(SCAL) -l $(INST)/sc$*/inst.list -d $(INST)/sc$*/ > $@

$(LOGS)/scalability.log: $(SCAL) $(SCAL_FILES)
	python $(SCAL) --header > $@ && \
	cat $(SCAL_FILES) >> $@ && \
	tar --remove-files -czf $(INST)/scalability_p$(SCAL_P)$(SCAL_R)x$(SCAL_K).tar.gz $(INST)/sc* && \
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
check_src:
	egrep -nr --color 'TODO|FIXME|BUG|NOTE'

# clean recipes
clean-raw-inst:
	@echo Cleaning up raw instances in $(INST)...
	rm -f $(INST)/R/*.bdd
	rm -f $(INST)/R/*.list*
	rm -f $(INST)/N/*.bdd
	rm -f $(INST)/N/*.list*
	rm -rf $(INST)/ds_stats
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
