# Aligning BDDs: the main makefile
# Contains all key instructions to reproduce the paper
# results (figures)
#
# (c) Alexey Bochkraev, Clemson University, 2020
# abochka@clemson.edu

######################################################################
## Files and directories
PREF=./
FIGS=./figures
INST=$(PREF)/instances
LOGS=$(PREF)/run_logs
PP=./post_processing

STATS=python -m gen_lsizes_stats
######################################################################
## Constants
JUFL_P=0.3
JUFL_N=20
TUFL_N=20
TUFL_P=0.3
RND_N=15
RND_P=0.6
HEU_BM_NO_INST=10

ORIG_N = 10
ORIG_K_TOTAL = 100

STAT_N = 10
STAT_PS = 0.2 0.5 0.8

SCAL_K=1

## Technical constants
# how many jobs to run in parallel?
PARFLAG=2

## Key targets
figures: \
	$(FIGS)/guessing.tex \
	$(FIGS)/simpl_heuristics.eps \
	$(FIGS)/LB.eps \
	$(FIGS)/orig_obj_histograms.eps \
	$(FIGS)/orig_runtimes.eps \
	$(FIGS)/.opts_stats \
	$(FIGS)/tUFLP_runtimes_overview.eps \
	$(FIGS)/various_simpl_vs_min.eps \
	$(FIGS)/orig_lwidth_stats.eps \
	$(FIGS)/tUFLP_runtimes_breakdown.eps

######################################################################
## Calculated values
HEU_BM_K=$(shell expr $(HEU_BM_NO_INST) / $(PARFLAG) + 1)
STAT_K=$(shell expr $(STAT_N) / $(PARFLAG) + 1)
ORIG_K=$(shell expr $(ORIG_K_TOTAL) / $(PARFLAG) + 1)


######################################################################
## Figures description

## The table (simplified problem heuristics)

$(FIGS)/guessing.tex: $(LOGS)/main_rnd_run.csv $(PP)/tab_guessing.R
				Rscript $(PP)/tab_guessing.R -i $< -o $@


## Comparison of the heuristics to the *simplified* problem

$(FIGS)/simpl_heuristics.eps: $(LOGS)/main_rnd_run.csv $(PP)/fig_simpl_heuristics.R
				Rscript $(PP)/fig_simpl_heuristics.R -i $< -o $@


## Different approaches to lower bound
$(FIGS)/LB.eps: $(LOGS)/simpl_LB.csv $(PP)/fig_LBs.R
				Rscript $(PP)/fig_LBs.R -i $< -o $@

$(LOGS)/simpl_LB.csv: $(INST)/orig_problem.list compare_simpl_LBs.py
				rm -f $<.tmp.* && \
				split -d -nl/$(PARFLAG) $< $<.LBs.tmp. && \
				parallel -j $(PARFLAG) python -m compare_simpl_LBs -l {} ">" $@.tmp.{#} -d $(INST)/orig_problem/ ::: $$(ls $<.LBs.tmp.*) && \
				python -m compare_simpl_LBs --header > $@ && \
				cat $@.tmp.* >> $@ && \
				rm $@.tmp.* && \
				rm $<.LBs.tmp.*


## Histograms of the original objective value
$(FIGS)/orig_obj_histograms.eps: $(LOGS)/main_rnd_run.csv $(PP)/fig_obj_hist.R
				Rscript $(PP)/fig_obj_hist.R -i $< -o $@

$(LOGS)/main_rnd_run.csv: $(INST)/orig_problem.list
				rm -f $<.tmp.* && \
				split -d -nl/$(PARFLAG) $< $<.tmp. && \
				parallel -j $(PARFLAG) python -m solve_inst -i {} -o $@.tmp.{#} -d $(INST)/orig_problem/ ::: $$(ls $<.tmp.*) && \
				python -m solve_inst --header > $@ && \
				cat $@.tmp.* >> $@ && \
				rm $@.tmp.* && \
				rm $<.tmp.*

$(INST)/orig_problem.list: gen_BDD_pair.py
				rm -rf $(INST)/orig_problem && \
				mkdir $(INST)/orig_problem && \
				parallel -j $(PARFLAG) 'python -m gen_BDD_pair -s $$(( $(ORIG_K) * {#} - $(ORIG_K)))' -K $(ORIG_K) -v $(ORIG_N) -p $(RND_P) -R -U $(INST)/orig_problem/ --quiet ::: $(shell seq 1 $(PARFLAG)) && \
				ls $(INST)/orig_problem | grep -Po "A\\K[^\\.]*" > $@

## Original problem (rnd) scaling figure
$(FIGS)/orig_runtimes.eps: $(LOGS)/orig_scal.csv $(PP)/fig_scal.R
				Rscript $(PP)/fig_scal.R -i $< -o $@


$(LOGS)/orig_scal.csv: $(INST)/scal/instances.list par_scal_test.py
				rm -f $(INST)/scal/instances.list.* && \
				split -d -nl/$(PARFLAG) $< $(INST)/scal/instances.list. && \
				parallel -j $(PARFLAG) python -m par_scal_test -i {} -o $@.tmp.{#} -d $(INST)/scal/ ::: $$(ls $(INST)/scal/instances.list.*) && \
				python -m par_scal_test --header > $@ && \
				cat $@.tmp.* >> $@ && \
				rm $@.tmp.* && \
				tar czf $(INST)/orig_scal.tar.gz -C $(INST)/scal . --remove-files


$(INST)/scal/instances.list: gen_BDD_pair.py
				mkdir -p $(INST)/scal && \
				parallel -j $(PARFLAG) 'python -m gen_BDD_pair -s $$(( $(SCAL_K) * {#} - $(SCAL_K)))' -K $(SCAL_K) -v {} -p $(RND_P) -R -U $(INST)/scal/ --quiet ::: $(shell seq 5 7 | shuf) && \
				ls $(INST)/scal | grep -Po "A\\K[^\\.]*" > $@


## Optima stats for the simplified problem
$(FIGS)/no_opts.eps $(FIGS)/opts_diam.eps $(FIGS)/heuristic_simscore.eps $(FIGS)/heuristic_simscore_vs_AB_simscore.eps: $(FIGS)/.opts_stats

$(FIGS)/.opts_stats: $(LOGS)/simpl_sol_struct.csv $(PP)/figs_simpl_opt_struct.R
				Rscript $(PP)/figs_simpl_opt_struct.R -i $< --outdir $(FIGS) && \
				touch $(FIGS)/.opts_stats

$(LOGS)/simpl_sol_struct.csv: experiments/heu_sol_struct.py
				python -m experiments.heu_sol_struct -k 10 -n 6 > $@

## t-UFLP runtimes overview (scaling)
$(FIGS)/tUFLP_runtimes_overview.eps: $(LOGS)/tUFLP_runtimes_scal.csv $(PP)/fig_tUFLP_runtimes_scal.R
				Rscript $(PP)/fig_tUFLP_runtimes_scal.R -i $< -o $@

$(LOGS)/tUFLP_runtimes_scal.csv: experiments/tUFLP_runtimes.py
				parallel -j $(PARFLAG) python -m experiments.tUFLP_runtimes -K 100 -n {} --with-gsifts ">" $@.tmp.{} ::: $(shell seq 5 10 | shuf) && \
				python -m experiments.tUFLP_runtimes -H > $@ && \
				cat $@.tmp.* >> $@ && rm $@.tmp.*

$(FIGS)/various_simpl_vs_min.eps: $(PP)/fig_simpl_vs_min.R $(LOGS)/heu_bm/jUFLP.csv $(LOGS)/heu_bm/tUFLP_nat.csv $(LOGS)/heu_bm/tUFLP_rnd.csv $(LOGS)/heu_bm/rnd_dia.csv
				Rscript post_processing/fig_simpl_vs_min.R --indir $(LOGS)/heu_bm --out $@

######################################################################
## Logfiles generation

$(LOGS)/heu_bm/jUFLP.csv: experiments/jUFL_hist_sizes.py
				mkdir -p $(LOGS)/heu_bm && \
				python -m experiments.jUFL_hist_sizes -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.jUFL_hist_sizes -n $(JUFL_N) -K $(HEU_BM_K) -p $(JUFL_P) -P {} -l $(INST)/jUFLP_inst.tmp.{} ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				cat $(INST)/jUFLP_inst.tmp* > $(INST)/jUFLP_instances.json && \
				rm $(INST)/jUFLP_inst.tmp*


$(LOGS)/heu_bm/tUFLP_nat.csv: experiments/tUFL_hist_sizes_control.py
				mkdir -p $(LOGS)/heu_bm && \
				python -m experiments.tUFL_hist_sizes_control -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.tUFL_hist_sizes_control -n $(TUFL_N) -K $(HEU_BM_K) -p $(TUFL_P) -P {} -l $(INST)/tUFLP_nat_inst.tmp.{} ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				cat $(INST)/tUFLP_nat_inst.tmp* > $(INST)/tUFLP_nat_instances.json && \
				rm $(INST)/tUFLP_nat_inst.tmp*


$(LOGS)/heu_bm/tUFLP_rnd.csv: experiments/tUFL_hist_sizes_control.py
				mkdir -p $(LOGS)/heu_bm && \
				python -m experiments.tUFL_hist_sizes_control -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.tUFL_hist_sizes_control -R -n $(TUFL_N) -K $(HEU_BM_K) -p $(TUFL_P) -P {} -l $(INST)/tUFLP_rnd_inst.tmp.{} ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				cat $(INST)/tUFLP_rnd_inst.tmp* > $(INST)/tUFLP_rnd_instances.json && \
				rm $(INST)/tUFLP_rnd_inst.tmp*


$(LOGS)/heu_bm/rnd_dia.csv: experiments/rnd_dia_hist_sizes_control.py
				mkdir -p $(LOGS)/heu_bm && \
				rm -rf $(INST)/bm_heu_inst && \
				mkdir -p $(INST)/bm_heu_inst && \
				python -m experiments.rnd_dia_hist_sizes_control -H > $@ && \
				parallel -j $(PARFLAG) python -m experiments.rnd_dia_hist_sizes_control -n $(RND_N) -K $(HEU_BM_K) -p $(RND_P) -P {} -l $(INST)/bm_heu_inst ">>" $@.tmp.{} ::: $(shell seq 1 $(PARFLAG)) && \
				cat $@.tmp.* >> $@ && rm $@.tmp.* && \
				tar czf $(INST)/bm_heu_random_diagrams.tar.gz -C $(INST)/bm_heu_inst . --remove-files


## Random dataset stats
$(FIGS)/orig_lwidth_stats.eps: $(LOGS)/lwidths.log $(PP)/fig_summary.R
				Rscript $(PP)/fig_summary.R -i $< -o $@

$(LOGS)/lwidths.log: gen_BDD_pair.py
				mkdir -p $(INST)/orig_stats && \
				parallel mkdir -p $(INST)/orig_stats/{} ::: $(STAT_PS) && \
				parallel -j $(PARFLAG) python -m gen_BDD_pair -v $(RND_N) -K $(STAT_K) -p $(RND_P) -R -U $(INST)/orig_stats/{1} --quiet ::: $(STAT_PS) ::: $(shell seq 1 $(PARFLAG)) && \
				parallel find $(INST)/orig_stats/{}/*.bdd ">" $(INST)/orig_stats/{}/inst.list ::: $(STAT_PS) && \
				$(STATS) -H $@ && \
				parallel $(STATS) -s {} $(INST)/orig_stats/{}/inst.list ">" $(LOGS)/lwidths.tmp.{} ::: $(STAT_PS) && \
				cat $(LOGS)/lwidths.tmp.* >> $@ && rm $(LOGS)/lwidths.tmp.* && \
				tar czf $(INST)/orig_stats.tar.gz -C $(INST)/orig_stats . --remove-files


## Breakdown of runtimes (CPP)

$(FIGS)/tUFLP_runtimes_breakdown.eps: $(LOGS)/tUFLP_runtimes.csv $(PP)/fig_tUFLP_runtimes_breakdown.R
				Rscript $(PP)/fig_tUFLP_runtimes_breakdown.R -i $< -o $@

$(LOGS)/tUFLP_runtimes.csv: experiments/tUFLP_runtimes.py
				python -m experiments.tUFLP_runtimes -H > $@ && \
				python -m experiments.tUFLP_runtimes -K 20 -n 20 -l $(INST)/tUFLP_steps_breakdown.json >> $@

## clean-up
save-orig-instances:
				tar czf $(INST)/orig_problem.tar.gz -C $(INST)/orig_problem . --remove-files
######################################################################
## Files and directories
# INST=$(PREF)/instances/raw
# ARC=./data_archive
# FIGS=./figures
# SOLVE=./solve_inst.py
# BBLOG=./log_bb.py
# SCAL=./scal_test.py
# STATS=./gen_lsizes_stats.py
# CMP_LBs=./compare_simpl_LBs.py
# PSCAL=./par_scal_test.py
# ######################################################################
# ## Numerical parameters
# PAR_SOL=8

# ### random dataset stats
# LW_n=10000
# LW_N=15
# LW_Ps=0.2 0.4 0.6 0.8

# ### dataset generation parameters:
# p=0.6# dataset generation parameter
# n=10000# number of instances
# N=15# number of variables per instance
# n_LBs=1000

# ### scalability figure
# SCAL_N=5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 #26 27 28 29 30
# SCAL_P=$(p)
# SCAL_R=R

# ######################################################################
# ## calculated vars/parameters
# SCAL_FILES=$(addsuffix .log, $(addprefix $(LOGS)/scal_$(SCAL_R),$(SCAL_N)))
# LW_FILES_R=$(addsuffix .log, $(addprefix $(LOGS)/lwidths_Rp,$(LW_Ps)))
# LW_FILES_N=$(addsuffix .log, $(addprefix $(LOGS)/lwidths_Np,$(LW_Ps)))

# MAX_I=$(shell expr $(PAR_SOL) - 1 )
# SOL_FILES_R = $(addsuffix .log, $(addprefix $(LOGS)/part.solved_R., $(shell seq -f %02g 0 $(MAX_I))))
# SOL_FILES_N = $(addsuffix .log, $(addprefix $(LOGS)/part.solved_N., $(shell seq -f %02g 0 $(MAX_I))))

# BB_FILES_R = $(addsuffix .log, $(addprefix $(LOGS)/part.BB_bounds_R., $(shell seq -f %02g 0 $(MAX_I))))
# BB_FILES_N = $(addsuffix .log, $(addprefix $(LOGS)/part.BB_bounds_N., $(shell seq -f %02g 0 $(MAX_I))))

# PSCAL_FILES = $(addsuffix .log, $(addprefix $(LOGS)/part.scal., $(shell seq -f %02g 0 $(MAX_I))))

# DTE=$(shell date +%F)
# TS=$(shell date +%H-%M-%S)

# MAX_I=$(shell expr $(PAR_SOL) - 1 )

# .PHONY: all figures clean-raw-inst clean-insts move-logs figures/sample_BB_tree.png

# ######################################################################
# ## High-level recipes
# .SECONDARY: # keep all the intermediary logfiles (will not work o/w)

# all:

# figures/sample_BB_tree.png: ./sample_BB_tree.py
# 	python ./sample_BB_tree.py -v -V 8 -n 10 -o ./run_logs/sample_BB_tree.dot && \
# 	dot -Tpng ./run_logs/sample_BB_tree.dot > ./figures/sample_BB_tree.png


# ######################################################################
# ## Figure recipes
# SOL_RECIPE = Rscript $(PP)/fig_$*.R -i $< -o $@
# BB_RECIPE = Rscript $(PP)/fig_BB_$*.R -i $< -o $@


# $(FIGS)/fig_BB_%_R.eps: $(LOGS)/BB_bounds_R.log $(PP)/fig_BB_%.R
# 	$(BB_RECIPE)

# $(FIGS)/fig_BB_%_N.eps: $(LOGS)/BB_bounds_N.log $(PP)/fig_BB_%.R
# 	$(BB_RECIPE)


# ######################################################################
# ######################################################################


# $(LOGS)/BB_bounds_%.log: $$(BB_FILES_%)
# 	python $(BBLOG) --header > $@ && \
# 	tail -qn +2 $(BB_FILES_$*) >> $@

# $(LOGS)/part.BB_bounds_R.%.log: $(INST)/R/instances.list $(SOLVE)
# 	python $(BBLOG) -i $<.$* -o $@ -d $(INST)/R/

# $(LOGS)/part.BB_bounds_N.%.log: $(INST)/N/instances.list $(SOLVE)
# 	python $(BBLOG) -i $<.$* -o $@ -d $(INST)/N/


# ######################################################################
# # auxiliary recipes

prep_dirs:
				mkdir -p $(INST)
				mkdir -p $(LOGS)

save_data:
				cp -r $(INST)/* ./instances/
				cp -r $(LOGS)/* ./run_logs/

# check_src:
# 	egrep -nr --color 'TODO|FIXME|BUG|NOTE'

# install_R_pkg:
# 	Rscript ./aux/R_install_packages.R

# move-logs:
# 	@echo moving logs away from $(LOGS) to $(ARC)
# 	rm -f $(LOGS)/part.* && \
# 	tar --remove-files -czf $(ARC)/$(DTE)_$(TS)_logs.tar.gz -C $(LOGS) .

# # clean recipes
# clean-raw-inst:
# 	@echo Cleaning up raw instances in $(INST)...
# 	rm -f $(INST)/R/*.bdd
# 	rm -f $(INST)/R/*.list*
# 	rm -f $(INST)/N/*.bdd
# 	rm -f $(INST)/N/*.list*
# 	rm -rf $(INST)/ds_stats
# 	rm -f $(INST)/*.bdd

# clean-archive:
# 	@echo Cleaning up instance archives...
# 	rm -f $(ARC)/*.tar.gz

# clean-logs:
# 	@echo Cleaning up log files...
# 	rm -f $(LOGS)/*.log

# clean-figures:
# 	@echo Cleaning up figures...
# 	rm -f $(FIGS)/*.eps

# clean: clean-raw-inst clean-logs clean-figures
# clean-tmp: clean-raw-inst clean-logs

# test: $(LOGS)/BB.test
# 	@echo Running tests...

# $(LOGS)/BB.test: $(LOGS)/BB_bounds_R.log $(LOGS)/solved_R.log tests/BB_log_correct.R
# 	Rscript tests/BB_log_correct.R -b $(LOGS)/BB_bounds_R.log -s $(LOGS)/solved_R.log > $(LOGS)/BB.test

# qtest:
# 	python -m pytest tests/UFL_test.py
